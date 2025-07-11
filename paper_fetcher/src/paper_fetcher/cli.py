import argparse
import sys
from typing import List

from paper_fetcher.main import (
    search_pubmed,
    fetch_details,
    extract_paper_info,
    write_csv,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch PubMed papers with at least one industry-affiliated author."
    )
    parser.add_argument("query", type=str, help="PubMed query string (in quotes)")
    parser.add_argument(
        "-f", "--file", type=str, help="Output CSV file (if not provided, print to console)"
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enable debug logging"
    )

    args = parser.parse_args()

    if args.debug:
        print(f"[DEBUG] Searching PubMed for: {args.query}", file=sys.stderr)

    pmids = search_pubmed(args.query)
    # pmids = [
    #     # "40636562","40636355","40636320","40636313","40636157","40636108","40635758","40635452","40635373""40634661","40634456","40634170","40634102","40634064","40633919","40633571","40633204","40633069","40632851","40632736","40632654","40632557","40632549","40632333"
    #          ]

    if args.debug:
        print(f"[DEBUG] Found {len(pmids)} articles", file=sys.stderr)

    articles = fetch_details(pmids)
    results = []

    # print(articles)
    
    for article in articles:
        record = extract_paper_info(article)
        if record:
            results.append(record)

    if args.debug:
        print(f"[DEBUG] {len(results)} articles with non-academic authors", file=sys.stderr)

    write_csv(results, args.file)
