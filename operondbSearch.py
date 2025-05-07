from playwright.sync_api import sync_playwright
import pandas as pd
from bs4 import BeautifulSoup

def fetch_operon_data(query: str, include_all_rows: bool) -> pd.DataFrame:
    """
    Выполняет поиск оперонов на OperonDB для заданного запроса с использованием Playwright.

    Аргументы:
      query: Строка поиска, например "Haemophilus influenzae Rd KW20".
      include_all_rows:
          True  – сохранить все найденные строки;
          False – сохранять только строки с непустым Definition.
                  Если встречается первая строка с пустым Definition, обход прекращается.

    Возвращает:
      DataFrame с колонками: "Operon ID", "Species", "Name", "Genes" (список строк) и "Definition".
    """
    base_url = "https://operondb.jp/search"
    key = "species"
    page = 1
    data = []  # Здесь будем аккумулировать найденные строки таблицы
    header = None  # Заголовки столбцов
    formatted_query = query.replace(" ", "%20")

    with sync_playwright() as p:
        # Запускаем Chromium в headless-режиме
        browser = p.chromium.launch(headless=True)
        context = browser.new_context()
        p_page = context.new_page()

        while True:
            url = f"{base_url}?p={page}&q={formatted_query}&key={key}"
            print(f"Обработка страницы {page}: {url}")
            try:
                # Переходим по URL с таймаутом 30 сек.
                p_page.goto(url, timeout=30000)
                # Сначала ждем появления основного контейнера #app (примерно 15 сек.)
                p_page.wait_for_selector("#app", timeout=15000)
                # Даем дополнительное время (например, 5 сек.) для того, чтобы динамический контент загрузился.
                p_page.wait_for_timeout(5000)
                html = p_page.content()
            except Exception as e:
                print(f"Ошибка загрузки страницы {page}: {e}")
                break

            # Выводим длину полученного HTML для отладки
            print(f"Длина HTML (страница {page}): {len(html)}")

            soup = BeautifulSoup(html, "html.parser")

            # На первой странице извлекаем заголовки из элемента <thead> внутри #app, если он есть
            if page == 1:
                thead = soup.select_one("div#app thead") or soup.select_one("thead")
                if thead:
                    header = [th.get_text(strip=True) for th in thead.find_all("th")]
                if not header or len(header) < 5:
                    header = ["Operon ID", "Species", "Name", "Genes", "Definition"]

            # Пытаемся найти таблицу, которая находится внутри контейнера #app
            table = soup.select_one("div#app table")
            if not table:
                print(f"Таблица не найдена на странице {page}. Завершение обхода.")
                break

            rows = table.find_all("tr")
            # Если таблица содержит только строку заголовков или пуста – завершаем
            if not rows or len(rows) <= 1:
                print(f"Нет строк с данными в таблице на странице {page}. Завершение обхода.")
                break

            # Пропускаем первую строку (заголовки) и обрабатываем оставшиеся строки
            data_rows = rows[1:]
            new_data_count = 0
            stop_flag = False

            for row in data_rows:
                cells = row.find_all("td")
                if len(cells) < 5:
                    continue  # Пропускаем строки, у которых меньше пяти ячеек

                # Извлекаем данные по столбцам
                operon_id = cells[0].get_text(strip=True)
                species = cells[1].get_text(strip=True)
                name = cells[2].get_text(strip=True)
                genes_cell = cells[3]
                # Из столбца "Genes" собираем список текста из каждого <li>
                genes = [li.get_text(strip=True) for li in genes_cell.find_all("li")]
                definition = cells[4].get_text(strip=True)

                if not include_all_rows and definition == "":
                    print("Найдена строка с пустым Definition. Прерывание обхода.")
                    stop_flag = True
                    break

                data.append([operon_id, species, name, genes, definition])
                new_data_count += 1

            if new_data_count == 0:
                print(f"На странице {page} не найдено новых данных. Завершение обхода.")
                break

            if stop_flag:
                break

            page += 1

        browser.close()

    df = pd.DataFrame(data, columns=header)
    return df