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
    """
    Entry point for the command-line tool.

    Parses arguments, performs a PubMed search, filters for non-academic authors
    affiliated with pharmaceutical/biotech companies, and saves the results as CSV.
    """
    parser = argparse.ArgumentParser(
        description="Fetch PubMed papers with at least one industry-affiliated author."
    )

    # Required query argument (PubMed query string)
    parser.add_argument("query", type=str, help="PubMed query string (in quotes)")

    # Optional output CSV file
    parser.add_argument(
        "-f", "--file", type=str,
        help="Output CSV file (if not provided, print to console)"
    )

    # Optional debug flag to print extra information
    parser.add_argument(
        "-d", "--debug", action="store_true",
        help="Enable debug logging"
    )

    # Optional limit on number of articles (default 1000, max 9999)
    parser.add_argument(
        "-l", "--limit", type=int, default=1000,
        help="Maximum number of articles to fetch (1 to 9999, default: 1000)"
    )

    args = parser.parse_args()

    # Enforce limit boundaries (1–9999)
    if args.limit < 1 or args.limit > 9999:
        print("[ERROR] --limit must be between 1 and 9999", file=sys.stderr)
        sys.exit(1)

    if args.debug:
        print(f"[DEBUG] Query: {args.query}", file=sys.stderr)
        print(f"[DEBUG] Fetching up to {args.limit} articles", file=sys.stderr)

    try:
        # Step 1: Search PubMed and get PMIDs
        pmids = search_pubmed(args.query, max_results=args.limit)

        if not pmids:
            print("[WARNING] No articles found for the given query.", file=sys.stderr)
            sys.exit(0)

        if args.debug:
            print(f"[DEBUG] Found {len(pmids)} articles", file=sys.stderr)

        # Step 2: Fetch full article details
        articles = fetch_details(pmids)

        if not articles:
            print("[ERROR] Failed to retrieve article metadata from PubMed.", file=sys.stderr)
            sys.exit(1)

        # Step 3: Extract relevant fields from each article
        results = [extract_paper_info(a) for a in articles if extract_paper_info(a)]

        if args.debug:
            print(f"[DEBUG] {len(results)} articles with non-academic authors", file=sys.stderr)

        # Step 4: Write results to file or stdout
        write_csv(results, args.file)

    except Exception as e:
        print(f"[ERROR] An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)
