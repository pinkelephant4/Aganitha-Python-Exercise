from typing import List, Optional, Dict
from Bio import Entrez
import csv
import re
import sys
import time


Entrez.email = "ananya.219311187@muj.manipal.edu"
MAX_PUBMED_RESULTS = 9999


def search_pubmed(query: str, max_results: int = 1000) -> List[str]:
    pmids = []
    batch_size = 1000
    start = 0

    max_results = min(max_results, MAX_PUBMED_RESULTS)

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
        pmids.extend(batch_pmids)

        if not batch_pmids:
            break

        start += len(batch_pmids)
        time.sleep(0.34)

    return pmids


def fetch_details(pmids: List[str]) -> List[Dict]:
    records = []
    batch_size = 200

    for i in range(0, len(pmids), batch_size):
        batch_ids = pmids[i:i + batch_size]
        ids_str = ",".join(batch_ids)

        with Entrez.efetch(db="pubmed", id=ids_str, retmode="xml") as handle:
            fetched = Entrez.read(handle)

        records.extend(fetched.get("PubmedArticle", []))
        time.sleep(0.34)

    return records


def extract_paper_info(article: Dict) -> Optional[Dict]:

    # print(article)
    try:
        medline = article.get("MedlineCitation", {})
        # print(medline)
        # pmid = medline.get("PMID", {}).get("value", "")
        pmid = medline.get("PMID", {})
        # print(pmid)
        article_data = medline.get("Article", {})
        # print(article_data)
        title = article_data.get("ArticleTitle", "").strip()

        # print(pmid)
        # print(title)

        journal = article_data.get("Journal", {})
        pub_date = journal.get("JournalIssue", {}).get("PubDate", {})
        date_parts = [pub_date.get("Year", ""), pub_date.get("Month", ""), pub_date.get("Day", "")]
        publication_date = "-".join(part for part in date_parts if part)

        authors = article_data.get("AuthorList", [])
        non_academic_authors = []
        company_affiliations = []
        corresponding_email = ""

        company_keywords = ["Inc", "Ltd", "LLC", "Corp", "Corporation", "Co.", "Limited"]
        academic_keywords = ["University", "College", "Hospital", "Academy", "Unive", "Provincial", "National", "Nacional", "School", "Public"]

        for author in authors:
            aff_info = author.get("AffiliationInfo", [])
            if not aff_info:
                continue
            
            affiliation = aff_info[0].get("Affiliation", "")
            
            affiliation_lower = affiliation.lower()
            if not corresponding_email:
                match = re.search(r"[\w\.-]+@[\w\.-]+", affiliation)
                if match:
                    corresponding_email = match.group(0)

            has_company = any(keyword.lower() in affiliation_lower for keyword in company_keywords)
            has_academic = any(keyword.lower() in affiliation_lower for keyword in academic_keywords)

            if has_company and not has_academic:
                full_name = f"{author.get('ForeName', '').strip()} {author.get('LastName', '').strip()}"
                non_academic_authors.append(full_name)
                company_name = affiliation
                company_affiliations.append(company_name)

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
    except Exception:
        return None


def write_csv(records: List[Dict], filename: Optional[str] = None) -> None:
    fieldnames = [
        "PubmedID",
        "Title",
        "Publication Date",
        "Non-academic Authors",
        "Company Affiliations",
        "Corresponding Email"
    ]

    output_stream = open(filename, "w", newline="", encoding="utf-8") if filename else sys.stdout
    writer = csv.DictWriter(output_stream, fieldnames=fieldnames)
    writer.writeheader()

    for record in records:
        writer.writerow(record)

    if filename:
        output_stream.close()
