from playwright.sync_api import sync_playwright
import pandas as pd
from bs4 import BeautifulSoup

def fetch_operon_data(query: str, include_all_rows: bool) -> pd.DataFrame:
    """
    Searches for operons on OperonDB for the given query using Playwright.

    Arguments:
      query: Search string, for example "Haemophilus influenzae Rd KW20".
      include_all_rows:
          True  – save all found rows;
          False – save only rows with a non-empty Definition.
                  If the first row with an empty Definition is encountered, the search stops.

    Returns:
      DataFrame with columns: "Operon ID", "Species", "Name", "Genes" (list of strings) and "Definition".
    """
    base_url = "https://operondb.jp/search"
    key = "species"
    page = 1
    data = []  # Here we will accumulate the found table rows
    header = None  # Column headers
    formatted_query = query.replace(" ", "%20")

    with sync_playwright() as p:
        # Launch Chromium in headless mode
        browser = p.chromium.launch(headless=True)
        context = browser.new_context()
        p_page = context.new_page()

        while True:
            url = f"{base_url}?p={page}&q={formatted_query}&key={key}"
            print(f"Processing page {page}: {url}")
            try:
                # Navigate to the URL with a timeout of 30 seconds.
                p_page.goto(url, timeout=30000)
                # First, wait for the main container #app to appear (about 15 seconds)
                p_page.wait_for_selector("#app", timeout=15000)
                # Provide additional time (e.g., 5 seconds) for the dynamic content to load.
                p_page.wait_for_timeout(5000)
                html = p_page.content()
            except Exception as e:
                print(f"Error loading page {page}: {e}")
                break

            # Output the length of the retrieved HTML for debugging
            print(f"HTML length (page {page}): {len(html)}")

            soup = BeautifulSoup(html, "html.parser")

            # On the first page, extract headers from the <thead> element within #app, if it exists
            if page == 1:
                thead = soup.select_one("div#app thead") or soup.select_one("thead")
                if thead:
                    header = [th.get_text(strip=True) for th in thead.find_all("th")]
                if not header or len(header) < 5:
                    header = ["Operon ID", "Species", "Name", "Genes", "Definition"]

            # Try to find the table located inside the #app container
            table = soup.select_one("div#app table")
            if not table:
                print(f"Table not found on page {page}. Exiting loop.")
                break

            rows = table.find_all("tr")
            # If the table contains only the header row or is empty – exit the loop
            if not rows or len(rows) <= 1:
                print(f"No data rows found in table on page {page}. Exiting loop.")
                break

            # Skip the first row (headers) and process the remaining rows
            data_rows = rows[1:]
            new_data_count = 0
            stop_flag = False

            for row in data_rows:
                cells = row.find_all("td")
                if len(cells) < 5:
                    continue  # Skip rows that have fewer than five cells

                # Extract data from the columns
                operon_id = cells[0].get_text(strip=True)
                species = cells[1].get_text(strip=True)
                name = cells[2].get_text(strip=True)
                genes_cell = cells[3]
                # From the "Genes" column, collect a list of text from each <li>
                genes = [li.get_text(strip=True) for li in genes_cell.find_all("li")]
                definition = cells[4].get_text(strip=True)

                if not include_all_rows and definition == "":
                    print("Found a row with an empty Definition. Exiting loop.")
                    stop_flag = True
                    break

                data.append([operon_id, species, name, genes, definition])
                new_data_count += 1

            if new_data_count == 0:
                print(f"No new data found on page {page}. Exiting loop.")
                break

            if stop_flag:
                break

            page += 1

        browser.close()

    df = pd.DataFrame(data, columns=header)
    return df