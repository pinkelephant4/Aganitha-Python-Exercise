from typing import List, Optional, Dict
from Bio import Entrez
import csv
import re
import sys
import time

# Required by NCBI's Entrez API
Entrez.email = "ananya.219311187@muj.manipal.edu"  

# PubMed has a hard cap of 10,000 results per query
MAX_PUBMED_RESULTS = 9999


def search_pubmed(query: str, max_results: int = 1000) -> List[str]:
    """
    Searches PubMed for a given query and returns up to `max_results` PubMed IDs.

    Args:
        query (str): The PubMed query string.
        max_results (int): The maximum number of articles to retrieve (default 1000).

    Returns:
        List[str]: A list of PubMed IDs (PMIDs).
    """
    pmids = []
    batch_size = 1000
    start = 0

    max_results = min(max_results, MAX_PUBMED_RESULTS)

    try:
        while start < max_results:
            with Entrez.esearch(
                db="pubmed",
                term=query,
                retmode="xml",
                retstart=start,
                retmax=min(batch_size, max_results - start)
            ) as handle:
                results = Entrez.read(handle)

            batch_pmids = results.get("IdList", [])
            if not batch_pmids:
                break  # No more results

            pmids.extend(batch_pmids)
            start += len(batch_pmids)
            time.sleep(0.34)  # To respect NCBI rate limits

        return pmids

    except Exception as e:
        print(f"[ERROR] Failed to search PubMed: {e}", file=sys.stderr)
        return []


def fetch_details(pmids: List[str]) -> List[Dict]:
    """
    Fetches detailed article metadata from PubMed using a list of PMIDs.

    Args:
        pmids (List[str]): List of PubMed IDs to fetch metadata for.

    Returns:
        List[Dict]: A list of PubMed article records (raw).
    """
    records = []
    batch_size = 200

    for i in range(0, len(pmids), batch_size):
        batch_ids = pmids[i:i + batch_size]
        ids_str = ",".join(batch_ids)

        try:
            with Entrez.efetch(db="pubmed", id=ids_str, retmode="xml") as handle:
                fetched = Entrez.read(handle)

            records.extend(fetched.get("PubmedArticle", []))
            time.sleep(0.34)
        except Exception as e:
            print(f"[ERROR] Failed to fetch article details: {e}", file=sys.stderr)

    return records


def extract_paper_info(article: Dict) -> Optional[Dict]:
    """
    Extracts structured data from a PubMed article dictionary.

    Filters authors affiliated with companies and not academic institutions.

    Args:
        article (Dict): A single PubMed article record.

    Returns:
        Optional[Dict]: A dictionary of cleaned metadata or None if not relevant.
    """
    try:
        medline = article.get("MedlineCitation", {})
        article_data = medline.get("Article", {})
        pmid = medline.get("PMID", {})
        title = article_data.get("ArticleTitle", "").strip()

        # Extract publication date components
        journal = article_data.get("Journal", {})
        pub_date = journal.get("JournalIssue", {}).get("PubDate", {})
        date_parts = [pub_date.get("Year", ""), pub_date.get("Month", ""), pub_date.get("Day", "")]
        publication_date = "-".join(part for part in date_parts if part)

        authors = article_data.get("AuthorList", [])
        non_academic_authors = []
        company_affiliations = []
        corresponding_email = ""

        # Heuristic keywords
        company_keywords = ["Inc", "Ltd", "LLC", "Corp", "Corporation", "Co.", "Limited", "Pharma", "Biotech"]
        academic_keywords = [
            "University", "College", "Hospital", "Academy", "Unive", "Provincial",
            "National", "Nacional", "School", "Public"
        ]

        for author in authors:
            aff_info = author.get("AffiliationInfo", [])
            if not aff_info:
                continue

            affiliation = aff_info[0].get("Affiliation", "")
            affiliation_lower = affiliation.lower()

            # Try to extract an email from the affiliation
            if not corresponding_email:
                match = re.search(r"[\w\.-]+@[\w\.-]+", affiliation)
                if match:
                    corresponding_email = match.group(0)

            has_company = any(k.lower() in affiliation_lower for k in company_keywords)
            has_academic = any(k.lower() in affiliation_lower for k in academic_keywords)

            # Must be company-affiliated AND NOT academic-affiliated
            if has_company and not has_academic:
                full_name = f"{author.get('ForeName', '').strip()} {author.get('LastName', '').strip()}"
                non_academic_authors.append(full_name)
                company_affiliations.append(affiliation)

        if non_academic_authors:
            return {
                "PubmedID": pmid,
                "Title": title,
                "Publication Date": publication_date,
                "Non-academic Authors": "; ".join(non_academic_authors),
                "Company Affiliations": "; ".join(dict.fromkeys(company_affiliations)),
                "Corresponding Email": corresponding_email,
            }

        return None

    except Exception as e:
        print(f"[ERROR] Failed to extract data from article: {e}", file=sys.stderr)
        return None


def write_csv(records: List[Dict], filename: Optional[str] = None) -> None:
    """
    Writes a list of article records to a CSV file (or stdout if no file is given).

    Args:
        records (List[Dict]): List of extracted article data dictionaries.
        filename (Optional[str]): Output file path or None to print to stdout.
    """
    fieldnames = [
        "PubmedID",
        "Title",
        "Publication Date",
        "Non-academic Authors",
        "Company Affiliations",
        "Corresponding Email"
    ]

    try:
        output_stream = open(filename, "w", newline="", encoding="utf-8") if filename else sys.stdout
        writer = csv.DictWriter(output_stream, fieldnames=fieldnames)
        writer.writeheader()

        for record in records:
            writer.writerow(record)

        if filename:
            output_stream.close()

    except Exception as e:
        print(f"[ERROR] Failed to write CSV output: {e}", file=sys.stderr)
