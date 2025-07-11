from typing import List, Optional, Dict
from Bio import Entrez
import csv
import re
import sys

Entrez.email = "ananya.219311187@muj.manipal.edu"


def search_pubmed(query: str, retmax: int = 100) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    res = Entrez.read(handle)
    # print(res)
    return res.get("IdList", [])


def fetch_details(pmids: List[str]) -> List[Dict]:
    if not pmids:
        return []
    ids = ",".join(pmids)
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    res = Entrez.read(handle)
    # print(res)
    return res.get("PubmedArticle", [])


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
        academic_keywords = ["University", "College", "Hospital", "Academy", "Unive"]

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
                print(aff_info)
                print(affiliation)
                print()

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
